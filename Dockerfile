FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY dna_app ./dna_app
ENV DNA_DATA_DIR=/data
VOLUME /data
EXPOSE 8000
# Build the reference databases once with:
#   docker run --rm -v dna-data:/data IMAGE python -m dna_app.build_db clinvar
#   docker run --rm -v dna-data:/data IMAGE python -m dna_app.build_db 1000g
CMD ["uvicorn", "dna_app.server:app", "--host", "0.0.0.0", "--port", "8000"]
