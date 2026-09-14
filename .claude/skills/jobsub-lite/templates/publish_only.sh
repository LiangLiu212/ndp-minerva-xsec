#!/bin/bash
# Sentinel worker for the jobsub_lite publish flow.
# Runs as one real grid job so we can read $INPUT_TAR_FILE on the worker and
# report which /cvmfs/fifeuserN.opensciencegrid.org/sw/<group>/<hash>/ the RCDS
# publish actually landed in (N=1..4 is randomly assigned per upload).
echo "publish-only sentinel"
echo "PUBLISH_SENTINEL_PWD=$(pwd)"
echo "PUBLISH_SENTINEL_INPUT_TAR_FILE=${INPUT_TAR_FILE}"
echo "PUBLISH_SENTINEL_INPUT_TAR_DIR_LOCAL=${INPUT_TAR_DIR_LOCAL}"
if [ -n "${INPUT_TAR_FILE}" ]; then
  if [ -d "${INPUT_TAR_FILE}" ]; then
    echo "PUBLISH_SENTINEL_CVMFS_DIR=${INPUT_TAR_FILE%/}"
  else
    echo "PUBLISH_SENTINEL_CVMFS_DIR=$(dirname "${INPUT_TAR_FILE}")"
  fi
fi
exit 0
