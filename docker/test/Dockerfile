FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq -y && \
    apt-get install -qq -y  python3-pip \
                            python3-flask \
                            cifs-utils && \
    apt -y autoremove

ENV APP_USER="tracker" \
    APP_UID="1000" \
    APP_GROUP="tracker" \
    APP_GID="1000"

RUN set -ex; \
    groupadd -r --gid "$APP_GID" "$APP_GROUP"; \
    useradd -rm --uid "$APP_UID" --gid "$APP_GID" -d /home/"$APP_USER" "$APP_USER"

USER $APP_USER

RUN python3 -m venv --prompt asm-utils venv
RUN source venv/bin/activate
RUN pip install --upgrade pip
RUN pip install --editable .
ENV PATH=/home/$APP_USER/.local/bin:$PATH

