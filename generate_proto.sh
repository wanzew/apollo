#!/bin/bash 

files=$(find  ./ -name "*.proto")
set done_files=0
set total_files=0

for file in $files;do
  protoc -I . --cpp_out=. $file

  # files counter
  if [ $? == 0 ]; then
    echo $file "--done"
    done_files=$[$done_files+1]
    total_files=$[$total_files+1]
  elif [$? != 0 ]; then
    echo $file
    total_files=$[$total_files+1]
  fi
done

echo "total file:" $total_files
echo "done files:" $done_files