FROM uwgac/topmed-master:latest

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN apt-get update && apt-get -y install git

RUN git clone https://github.com/AnalysisCommons/LDGds.git