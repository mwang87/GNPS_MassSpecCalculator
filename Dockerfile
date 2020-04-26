FROM continuumio/miniconda3:latest
MAINTAINER Mingxun Wang "mwang87@gmail.com"

RUN apt-get update
RUN apt-get install -y build-essential
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . /app
WORKDIR /app

