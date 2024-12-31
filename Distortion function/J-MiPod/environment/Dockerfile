# hash:sha256:b3acb599635270d51eae193bbfe3c0376048369390d1ab6a99b7fba2068f1df7
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.4ubuntu1 \
        libjpeg62-dev=1:6b2-3 \
    && rm -rf /var/lib/apt/lists/*
