FROM michaelwetter/ubuntu-1604_jmodelica_trunk

USER root

MAINTAINER xxx

# This is not ideal. The tarballs are not named nicely and EnergyPlus versioning is strange
ENV ENERGYPLUS_VERSION 8.4.0
ENV ENERGYPLUS_TAG v8.4.0-Update1
ENV ENERGYPLUS_SHA 09f5359d8a

# This should be x.y.z, but EnergyPlus convention is x-y-z
ENV ENERGYPLUS_INSTALL_VERSION 8-4-0

# Downloading from Github
# e.g. https://github.com/NREL/EnergyPlus/releases/download/v8.4.0-Update1/EnergyPlus-8.4.0-09f5359d8a-Linux-x86_64.sh
ENV ENERGYPLUS_DOWNLOAD_BASE_URL https://github.com/NREL/EnergyPlus/releases/download/$ENERGYPLUS_TAG
ENV ENERGYPLUS_DOWNLOAD_FILENAME EnergyPlus-$ENERGYPLUS_VERSION-$ENERGYPLUS_SHA-Linux-x86_64.sh
ENV ENERGYPLUS_DOWNLOAD_URL $ENERGYPLUS_DOWNLOAD_BASE_URL/$ENERGYPLUS_DOWNLOAD_FILENAME


# Collapse the update of packages, download and installation into one command
# to make the container smaller & remove a bunch of the auxiliary apps/files
# that are not needed in the container

RUN apt-get update && apt-get install -y ca-certificates curl \
    && rm -rf /var/lib/apt/lists/* \
    && curl -SLO $ENERGYPLUS_DOWNLOAD_URL \
    && chmod +x $ENERGYPLUS_DOWNLOAD_FILENAME \
    && echo "y\r" | ./$ENERGYPLUS_DOWNLOAD_FILENAME \
    && rm $ENERGYPLUS_DOWNLOAD_FILENAME \
    && cd /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION \
    && rm -rf DataSets Documentation ExampleFiles WeatherData MacroDataSets PostProcess/convertESOMTRpgm \
    PostProcess/EP-Compare PreProcess/FMUParser PreProcess/ParametricPreProcessor PreProcess/IDFVersionUpdater

RUN pip install numpy

RUN pip install pandas

RUN pip install sklearn

# Remove the broken symlinks
RUN cd /usr/local/bin \
    && find -L . -type l -delete

# Add in the test files
#ADD test /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION/test_run
#RUN cp /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION/Energy+.idd \
#        /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION/test_run/

# Add a symbolink to Energy+.idd
RUN ["ln", "-s", "/usr/local/EnergyPlus-8-4-0/Energy+.idd", "/usr/local/Energy+.idd"]

VOLUME /var/simdata
WORKDIR /var/simdata


CMD [ "/bin/bash" ]
