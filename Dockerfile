FROM python:3.8

RUN apt-get update
RUN git clone https://github.com/mugundhan1/iLiSA.git && \
	cd iLiSA && \
	pip install .[dev] && \
	cd ~

RUN python -c "import ilisa; print(ilisa.__version__)"

