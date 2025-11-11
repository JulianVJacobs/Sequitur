#!/bin/sh

chmod 600 ${SSH_DIR}/config
chmod 600 ${SSH_DIR}/config.d/github
chmod 600 ${SSH_DIR}/github
chmod 644 ${SSH_DIR}/github.pub

eval "$(ssh-agent -s)"
ssh-add ${SSH_DIR}/github