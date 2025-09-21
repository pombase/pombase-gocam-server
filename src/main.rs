// A mini server for GO-CAMs

use getopts;

use axum::{
    extract::{Path, Request, State}, http::{header, StatusCode}, response::{IntoResponse, Response}, routing::get, Json, Router, ServiceExt
};

use pombase_gocam_process::{find_holes, model_connections_to_cytoscope, model_to_cytoscape_simple, GoCamCytoscapeStyle};
use serde_json::{json, Value};

use tracing_subscriber::EnvFilter;
use tokio::fs::read;
use tokio::io::AsyncReadExt as _;

use tower::{layer::Layer, timeout::TimeoutLayer};

use tower_http::normalize_path::NormalizePathLayer;
use tower_http::trace::TraceLayer;

use std::{collections::{BTreeSet, HashMap, HashSet}, io::Cursor, fs, process, sync::Arc, time::Duration};

use getopts::Options;

use pombase_gocam::{parse_gocam_model, GoCamGeneIdentifier, GoCamModel, GoCamModelId, GoCamNode, GoCamNodeOverlap, RemoveType};

type GoCamModelMap = HashMap<GoCamModelId, GoCamModel>;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} [options]", program);
    eprint!("{}", opts.usage(&brief));
}

struct AllState {
    web_root_dir: String,
    model_dir: String,
    gocam_models_by_id: GoCamModelMap,
    overlaps: Vec<GoCamNodeOverlap>,
    missing_activities: Vec<GoCamNode>,
}

fn read_gocam_models_from_dir(model_dir: &str)
    -> anyhow::Result<Vec<GoCamModel>>
{
    let mut ret = vec![];

    let entries = fs::read_dir(model_dir)?;

    for entry in entries {
        let file_name = entry?.path();
        if file_name.to_string_lossy().ends_with(".json") {
            let mut file = fs::File::open(file_name)?;
            let model = parse_gocam_model(&mut file)?;
            ret.push(model);
        }
    }

    Ok(ret)
}

async fn read_gocam_model(model_dir: &str, gocam_id: &str, flags: &HashSet<String>)
    -> anyhow::Result<GoCamModel>
{
    let file_name = format!("{}/{}.json", model_dir, gocam_id);
    let mut source = tokio::fs::File::open(file_name).await?;
    let mut contents = vec![];
    source.read_to_end(&mut contents).await?;
    let mut cursor = Cursor::new(contents);
    let model_res = parse_gocam_model(&mut cursor);

    let mut remove_types = HashSet::new();

    if flags.contains("no_inputs") {
        remove_types.insert(RemoveType::Targets);
    }
    if flags.contains("no_chemicals") {
        remove_types.insert(RemoveType::Chemicals);
    }

    if remove_types.is_empty() {
        model_res
    } else {
        model_res.map(|model| model.remove_nodes(remove_types))
    }
}

async fn read_gocam_models_by_ids(model_dir: &str,
                                  gocam_models_by_id: &GoCamModelMap)
    -> anyhow::Result<Vec<GoCamModel>>
{
    let mut models = vec![];

    let mut flags = HashSet::new();

    flags.insert("with_chemicals".to_owned());
    flags.insert("with_inputs".to_owned());

    for gocam_id in gocam_models_by_id.keys() {
        let model = read_gocam_model(model_dir, gocam_id, &flags).await?;
        models.push(model);
    }

    Ok(models)
}

pub async fn read_connected_gocam_models(model_dir: &str,
                                         gocam_models_by_id: &GoCamModelMap,
                                         overlaps: &Vec<GoCamNodeOverlap>,
                                         flags: &HashSet<String>)
    -> anyhow::Result<GoCamModel>
{
    if gocam_models_by_id.is_empty() {
        panic!();
    }

    let mut overlapping_gocam_ids = HashSet::new();
    for overlap in overlaps {
        for (model_id, _, _) in overlap.models.iter() {
            overlapping_gocam_ids.insert(model_id.replace("gomodel:", ""));
        }
    }

    let mut models = vec![];

    for gocam_id in gocam_models_by_id.keys() {
        if !overlapping_gocam_ids.contains(gocam_id.as_str()) {
            continue;
        }

        let flags = HashSet::new();
        let model = read_gocam_model(model_dir, gocam_id, &flags).await?;
        models.push(model);
    }

    let merge_res = GoCamModel::merge_models("merged", "merged models", &models);
    let mut remove_types = HashSet::new();

    if flags.contains("no_chemicals") {
        remove_types.insert(RemoveType::Chemicals);
    }
    if flags.contains("no_inputs") {
        remove_types.insert(RemoveType::Targets);
    }

    if remove_types.is_empty() {
        merge_res
    } else {
        merge_res.map(|model| model.remove_nodes(remove_types)
                         .retain_largest_subgraph())
    }
}

async fn read_merged_gocam_model(model_dir: &str,
                                 gocam_models_by_id: &GoCamModelMap,
                                 flags: &HashSet<String>,
                                 gene_list: &BTreeSet<GoCamGeneIdentifier>)
   -> anyhow::Result<GoCamModel>
{
    if gocam_models_by_id.is_empty() {
        panic!();
    }

    let models: Vec<_> = read_gocam_models_by_ids(model_dir, gocam_models_by_id).await?
        .into_iter()
        .filter(|model| {
            if flags.contains("trim_models") {
                model.model_activity_enabled_by(gene_list)
            } else {
                true
            }
        })
        .collect();

    let merge_res = GoCamModel::merge_models("merged", "merged models", &models);

    let mut remove_types = HashSet::new();

    if flags.contains("no_chemicals") {
        remove_types.insert(RemoveType::Chemicals);
    }
    if flags.contains("no_inputs") {
        remove_types.insert(RemoveType::Targets);
    }

    if remove_types.is_empty() {
        merge_res
    } else {
        merge_res.map(|model| model.remove_nodes(remove_types))
    }
}

async fn not_found() -> Json<Value> {
    json!({
        "status": "error",
        "reason": "Resource was not found."
    }).into()
}

async fn get_static_file(path: &str) -> Response {
    let res = read(path).await;

    let content_type = match mime_guess::from_path(&path).first_raw() {
        Some(mime) => mime,
        None => "text/plain"
    };

    match res {
        Ok(bytes) => {
            (StatusCode::OK, [(header::CONTENT_TYPE, content_type.to_string())], bytes).into_response()
        },
        Err(_) => {
            (StatusCode::NOT_FOUND, [(header::CONTENT_TYPE, "text/plain".to_string())], "not found".to_string()).into_response()
        }
    }
}

async fn get_index(State(all_state): State<Arc<AllState>>) -> Response {
    let web_root_dir = &all_state.web_root_dir;
    get_static_file(&format!("{}/index.html", web_root_dir)).await
}

async fn get_overlaps(State(all_state): State<Arc<AllState>>)
     -> impl IntoResponse
{
    let overlaps = &all_state.overlaps;

    Json(overlaps.to_owned())
}

async fn get_missing_activities(State(all_state): State<Arc<AllState>>)
     -> impl IntoResponse
{
    let holes = &all_state.missing_activities;

    Json(holes.to_owned())
}


async fn get_cytoscape_gocam_by_id(Path(gocam_id_arg): Path<String>,
                                   State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    get_cytoscape_gocam_by_id_retain_genes(Path((gocam_id_arg, "".to_owned())), State(all_state)).await
}

async fn get_cytoscape_gocam_by_id_retain_genes(Path((gocam_id_arg, gene_list)): Path<(String, String)>,
                                                State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let gocam_models_by_id = &all_state.gocam_models_by_id;

    if gocam_models_by_id.is_empty() {
        panic!();
    }

    let overlaps = &all_state.overlaps;
    let model_dir = &all_state.model_dir;

    let id_split: Vec<_> = gocam_id_arg.split(':').collect();
    let (gocam_id, flag_string) = (id_split[0], id_split.get(1));
    let flags: HashSet<_> =
        if let Some(flag_string) = flag_string {
            flag_string.split(",").map(String::from).collect()
        } else {
            HashSet::new()
        };
    let gene_set: BTreeSet<_> =
        gene_list.split(",").map(|g| g.to_owned()).collect();

    let mut read_res =
        if gocam_id.contains("+") {
            let filtered_data = gocam_id.split("+")
               .filter_map(|gocam_id| {
                    if let Some(model) = gocam_models_by_id.get(gocam_id) {
                        Some((gocam_id.into(), model.clone()))
                    } else {
                        None
                    }
                }).collect();
            read_merged_gocam_model(model_dir, &filtered_data,
                                    &flags, &gene_set).await
        } else {
            match gocam_id {
                "ALL_MERGED" => {
                    read_merged_gocam_model(model_dir, gocam_models_by_id,
                                            &flags, &gene_set).await
                }
                "ALL_CONNECTED" => {
                    read_connected_gocam_models(model_dir, gocam_models_by_id,
                                                overlaps,
                                                &flags).await
                }
                _ => read_gocam_model(model_dir, &gocam_id,
                                      &flags).await
            }
        };

    let style = if flags.contains("hide_models") {
        GoCamCytoscapeStyle::HideParents
    } else {
        GoCamCytoscapeStyle::IncludeParents
    };

    if flags.contains("retain_genes") {
        // we're hightlighting genes in the mega-models then hiding
        // models that don't have a hightlighted gene
        read_res = read_res.map(|model| {
            let mut remove_types = HashSet::new();
            remove_types.insert(RemoveType::Chemicals);
            remove_types.insert(RemoveType::Targets);
            model.remove_nodes(remove_types)
                .retain_enabling_genes(&gene_set)
        });
    }

    match read_res {
        anyhow::Result::Ok(model) => {
            let elements =  model_to_cytoscape_simple(&model, overlaps, style);

         (StatusCode::OK, Json(elements)).into_response()
        },
        anyhow::Result::Err(err) => {
            eprintln!("err: {}", err);
           (StatusCode::NOT_FOUND, [(header::CONTENT_TYPE, "text/plain".to_string())], "not found".to_string()).into_response()
        }
    }
}

async fn get_model_summary_for_cytoscape_all(State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let overlaps = &all_state.overlaps;

    let gocam_models_by_id = &all_state.gocam_models_by_id;

    let ids_and_titles: Vec<(String, String)> =
        gocam_models_by_id.iter()
        .map(|(gocam_id, model)| {
            (format!("gomodel:{}", gocam_id), model.title().to_owned())
        })
        .collect();

    let model_connections = model_connections_to_cytoscope(overlaps, &ids_and_titles);

    Json(model_connections.to_owned())
}

async fn get_model_summary_for_cytoscape_connected(State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let overlaps = &all_state.overlaps;

    let model_connections = model_connections_to_cytoscope(overlaps, &vec![]);

    Json(model_connections.to_owned())
}


#[tokio::main]
async fn main() {
    println!("{} v{}", PKG_NAME, VERSION);

    let args: Vec<String> = std::env::args().collect();
    let mut opts = Options::new();

    opts.optflag("h", "help", "print this help message");
    opts.optopt("b", "bind-address-and-port", "The address:port to bind to", "BIND_ADDRESS_AND_PORT");
    opts.optopt("w", "web-root-dir", "Root web data directory", "WEB_ROOT_DIR");
    opts.optopt("g", "gocam-model-directory",
                "The directory containing the GO-CAM model JSON files", "DIR");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("Invalid options\n{}", f)
    };

    let program = args[0].clone();

    if matches.opt_present("help") {
        print_usage(&program, &opts);
        process::exit(0);
    }

    let web_root_dir = matches.opt_str("w").unwrap_or_else(|| {
        eprintln!("no --web-root-dir|-w option");
        print_usage(&program, &opts);
        process::exit(1);
    });

    let gocam_model_dir = matches.opt_str("gocam-model-directory").unwrap_or_else(|| {
        eprintln!("no --gocam-model-directory|-g option");
        print_usage(&program, &opts);
        process::exit(1);
    });

    let gocam_models = read_gocam_models_from_dir(&gocam_model_dir)
        .expect(&format!("failed to read gocam models from {}", gocam_model_dir));

    let gocam_models_by_id = gocam_models.iter()
        .map(|m| (m.id().strip_prefix("gomodel:").unwrap_or(m.id()).to_owned(), m.clone()))
        .collect();

    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env()
                .or_else(|_| EnvFilter::try_new("pombase-chado-json=warn,tower_http=warn"))
                .unwrap(),
        )
        .init();

    let bind_address_and_port = matches.opt_str("bind-address-and-port");
    let listener =
        if let Some(bind_address_and_port) = bind_address_and_port {
            tokio::net::TcpListener::bind(bind_address_and_port).await.unwrap()
        } else {
            tokio::net::TcpListener::bind("0.0.0.0:9500").await.unwrap()
        };

    let overlaps = GoCamModel::find_overlaps(&gocam_models);

    let missing_activities = gocam_models.iter()
        .flat_map(|m| find_holes(m)).collect();

    let all_state = AllState {
        web_root_dir,
        gocam_models_by_id,
        overlaps,
        missing_activities,
        model_dir: gocam_model_dir,
    };

    println!("Starting server ...");
    let app = Router::new()
        .route("/", get(get_index))
        .route("/cytoscape/model_summary/all_models", get(get_model_summary_for_cytoscape_all))
        .route("/cytoscape/model_summary/connected_only", get(get_model_summary_for_cytoscape_connected))
        .route("/cytoscape/model_view/{gocam_id}", get(get_cytoscape_gocam_by_id))
        .route("/cytoscape/model_view/{gocam_id}/{retain_genes}", get(get_cytoscape_gocam_by_id_retain_genes))
        .route("/missing_activities", get(get_missing_activities))
        .route("/overlaps", get(get_overlaps))
        .fallback(not_found)
        .with_state(Arc::new(all_state))
        .layer(TraceLayer::new_for_http());

    let app = NormalizePathLayer::trim_trailing_slash().layer(app);

    axum::serve(listener, ServiceExt::<Request>::into_make_service(app))
        .await
        .unwrap()
        .layer(TimeoutLayer::new(Duration::from_secs(120)));

    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::DEBUG)
        .init();
}
