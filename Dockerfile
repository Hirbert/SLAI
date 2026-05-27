#Use official Python base image
FROM swr.cn-north-4.myhuaweicloud.com/ddn-k8s/docker.io/python:3.9-slim

# Set working directory
WORKDIR /app

# Set environment variables (optional)
ENV PYTHONDONTWRITEBYTECODE=1 \
PYTHONUNBUFFERED=1

# Copy dependency file and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy all script files
COPY *.py /app/

# Use bash as entrypoint for convenient execution of different scripts
CMD ["/bin/bash"]
