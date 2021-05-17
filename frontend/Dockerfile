# Dockerfile

FROM python:3-buster
# FROM python:3
# ENV PYTHONUNBUFFERED 1

# install gninx and vim
RUN apt-get update && apt-get install nginx vim -y --no-install-recommends
COPY nginx.default /etc/nginx/sites-available/default
RUN ln -sf /dev/stdout /var/log/nginx/access.log \
    && ln -sf /dev/stderr /var/log/nginx/error.log

RUN mkdir -p /opt/app
RUN mkdir -p /opt/app/batterycycling
RUN mkdir -p /opt/app/gibbsml

COPY requirements.txt start-server.sh /opt/app/
COPY batterycycling /opt/app/batterycycling
COPY gibbsml /opt/app/gibbsml

WORKDIR /opt/app
RUN pip install -r requirements.txt
RUN pip install ./gibbsml
RUN chown -R www-data:www-data /opt/app

EXPOSE 8020
STOPSIGNAL SIGTERM
CMD ["/opt/app/start-server.sh"]
