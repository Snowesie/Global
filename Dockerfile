FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY dna_app ./dna_app
COPY start.sh .
ENV DNA_DATA_DIR=/data
EXPOSE 8000
CMD ["sh", "start.sh"]
