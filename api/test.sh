#!/bin/sh
curl -X GET "http://127.0.0.1:8000/qc?min.features=200&max.features=5000&max.mtpercent=5" -H "accept: text/plain"
