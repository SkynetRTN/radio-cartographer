#!/bin/bash

# Configuration
IMAGE_NAME="radio-cartographer"
CONTAINER_NAME="radio-cartographer-sandbox"
# Mapping the specific WORKDIR defined in your Dockerfile 
WORK_DIR="/skynet/radio-cartographer"

export MSYS2_ARG_CONV_EXCL="*"
export MSYS_NO_PATHCONV=1

# 1. BUILD CHECK
# Check if the image exists; if not, build it using the Dockerfile
if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
  echo "building image $IMAGE_NAME..."
  docker build -t $IMAGE_NAME .
fi

# 1.5. FRESHNESS CHECK
# Check if the running container is using the latest image.
if [ ! -z "$(docker ps -qf "name=^${CONTAINER_NAME}$")" ]; then
    CURRENT_IMAGE_ID=$(docker inspect --format '{{.Image}}' $CONTAINER_NAME)
    LATEST_IMAGE_ID=$(docker inspect --format '{{.Id}}' $IMAGE_NAME)

    if [ "$CURRENT_IMAGE_ID" != "$LATEST_IMAGE_ID" ]; then
        echo "⚠️  Container is running an old image. Restarting..."
        docker rm -f $CONTAINER_NAME
    fi
fi

# 2. RUN CHECK (The "Persistent Sandbox")
# Check if container is running. If not, start it in "Sleep Mode".
# We use 'tail -f /dev/null' to keep the C++ environment alive indefinitely.
if [ -z "$(docker ps -qf "name=^${CONTAINER_NAME}$")" ]; then
  echo "⚠️  Sandbox container not running. Starting..."

  # VOLUME MOUNT CRITICAL:
  # -v "$(pwd):$WORK_DIR" maps your host source code to the container.
  # This overwrites the 'COPY' instruction from the Dockerfile
  # allowing Antigravity to edit files locally while compiling remotely.
  docker run -d \
    --rm \
    --name "$CONTAINER_NAME" \
    -p 5678:5678 \
    -v "$(pwd):$WORK_DIR" \
    -v "$WORK_DIR/build" \
    -w "$WORK_DIR" \
    "$IMAGE_NAME" \
    tail -f /dev/null

  # Wait for container spin-up
  sleep 2
fi

# 3. TTY DETECTION (The Antigravity Fix)
# -t 1 checks if the input is a human terminal.
# If Human: use -it (Interactive + TTY) for colors and shell prompting.
# If Agent: use -i (Interactive only) to prevent "input device is not a TTY" errors.
if [ -t 1 ]; then
    DOCKER_FLAGS="-it"
else
    DOCKER_FLAGS="-i"
fi

# 4. EXECUTION
# If arguments are provided (e.g. 'make'), run them.
# If no arguments are provided, open a bash shell.
if [ $# -eq 0 ]; then
    docker exec $DOCKER_FLAGS "$CONTAINER_NAME" /bin/bash
else
    # Support "sh -c" style invocation by checking for -c
    if [ "$1" == "-c" ]; then
        shift
    fi
    docker exec $DOCKER_FLAGS "$CONTAINER_NAME" /bin/bash -c "$*"
fi
