#!/bin/sh
# Container entrypoint. With DNA_AUTO_BUILD=1, any missing reference database
# is downloaded and built in the background on boot (the UI shows progress and
# holds uploads until it finishes); the web server starts immediately either way.
set -e
if [ "$DNA_AUTO_BUILD" = "1" ]; then
  python -m dna_app.build_db missing &
fi
exec uvicorn dna_app.server:app --host 0.0.0.0 --port "${PORT:-8000}"
