#!/bin/bash 

files=$(find  ./ -name "*.proto")

for file in $files;do
  protoc -I . --cpp_out=. $file
  echo $file
done
