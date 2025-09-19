# build with:
# (cd ~/git/pombase-gocam-server
#  VERSION=`perl -ne 'print "$1" if /version *= *"(\d+\.\d+\.\d+)"$/' Cargo.toml`
#  docker build --build-arg version=$VERSION -f etc/Dockerfile -t=pombase/gocam-server:v$VERSION .)

FROM debian:trixie

LABEL maintainer="Kim Rutherford <kim@pombase.org>"

ARG version=unset

RUN echo "building v$version"

ENV mod_url=UNSET

RUN mkdir /data

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    apt-get -y install apt-utils curl git-core

RUN apt install -y python3-django python3-tornado python3-regex pipx gunicorn

RUN pipx install circus

RUN apt install -y jq rustc cargo

ENV foo=bar
RUN git clone https://github.com/pombase/pombase-python-web.git
RUN git clone https://github.com/pombase/pombase-gocam-server.git

RUN (curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain stable && \
     . "$HOME/.cargo/env" && \
     (cd pombase-gocam-server; rustup default stable; cargo build --profile releasedev) && \
     rm -r $HOME/.rustup)

HEALTHCHECK --start-period=20s --interval=10s --timeout=20s --retries=10 \
  CMD /pombase-gocam-server/etc/healthcheck

CMD /pombase-gocam-server/etc/start_all.sh

# server setup
# docker service create --publish 9999:9999 --publish 9400:9400 --env mod_url="http://flybase.org" --name gocam-server-fly --replicas 1 pombase/gocam-server:v0.3.0
