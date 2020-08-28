FROM python:3
ENV PYTHONUNBUFFERED 1

RUN mkdir /calc
WORKDIR /calc

COPY requirements.txt /calc/
RUN pip install -r requirements.txt

COPY ./batterycycling /calc/