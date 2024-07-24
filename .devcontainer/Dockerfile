FROM alpine

RUN apk fix && apk --no-cache --update add git build-base python3-dev py3-pip linux-headers 
RUN python3 -m venv /workspaces/Sequitur/.venv && source .venv/bin/activate && \
    pip install --upgrade pip && pip install -r /workspaces/Sequitur/.venv/requirements.txt