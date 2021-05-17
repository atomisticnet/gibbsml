# Dockerfile

FROM python:3-buster

RUN mkdir -p /opt/app
RUN mkdir -p /opt/app/batterycycling
RUN mkdir -p /opt/app/gibbsml

COPY requirements.txt start-server.sh /opt/app/
COPY batterycycling /opt/app/batterycycling
COPY gibbsml /opt/app/gibbsml

# install gninx and vim
RUN apt-get update && apt-get install nginx vim -y --no-install-recommends
RUN ln -sf /dev/stdout /var/log/nginx/access.log \
    && ln -sf /dev/stderr /var/log/nginx/error.log
COPY nginx.default /etc/nginx/sites-available/default

# install python packages
WORKDIR /opt/app
RUN pip install -r requirements.txt
RUN pip install ./gibbsml

# make sure permissions are set correctly
RUN chown -R www-data:www-data /opt/app

# expose port 8020
EXPOSE 8020

# quit when the SIGTERM signal is received (CTRL+C)
STOPSIGNAL SIGTERM

# run the start-server.sh script without shell
CMD ["/opt/app/start-server.sh"]
