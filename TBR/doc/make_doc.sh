#!/bin/bash
doxygen
cd ../doc/html
google-chrome index.html
