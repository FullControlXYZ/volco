FROM python:3.8.10

RUN pip3 install --upgrade pip

RUN mkdir /usr/src/volco
WORKDIR /usr/src/volco
COPY ./requirements.txt .
RUN pip3 install -r requirements.txt
ENV PYTHONUNBUFFERED 1
COPY . .
