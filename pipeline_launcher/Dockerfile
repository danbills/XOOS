FROM python:3.13-slim AS builder-base

RUN apt-get -qq update && apt-get -qq install -y --no-install-recommends \
      curl \
      git \
      git-lfs \
      procps \
      unzip \
  && rm -rf /var/lib/apt/lists/*

COPY LICENSE.md /LICENSE.md
