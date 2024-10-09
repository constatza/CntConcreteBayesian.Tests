#!/bin/bash

# Upload the files to the server
rsync -avzP ./bin/x64/Release/net8.0/linux-x64/ meluxina:apps/mat/bin/ 
rsync -avzP ./DataFiles/ meluxina:apps/mat/DataFiles/
