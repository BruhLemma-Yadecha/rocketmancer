FROM python:3

ENV PYTHONDONOTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

WORKDIR /app

COPY requirements/ /app/requirements/
RUN pip install --upgrade pip
RUN pip install -r requirements/development.txt

COPY . /app/
EXPOSE 8000
