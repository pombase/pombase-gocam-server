FROM debian:trixie

LABEL maintainer="Kim Rutherford <kim@pombase.org>"

ENV mod_url=UNSET

RUN mkdir /data

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    apt-get -y install apt-utils curl git-core

RUN apt install -y python3-django python3-django-allauth python3-tornado \
    pipx gunicorn python3-seaborn python3-numpy python3-pandas \
    python3-matplotlib python3-regex

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
