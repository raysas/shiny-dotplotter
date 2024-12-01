#!/bin/bash

pip install shiny
pip install rsconnect-python

rsconnect add \
	  --account raysas \
	  --name raysas \
	  --token AA55621F1C938E6862AD4DE670620779 \
	  --secret <SECRET> # replace with secret