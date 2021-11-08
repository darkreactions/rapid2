#!/usr/bin/env bash


# Usage:
# . compress_video.sh {FOLDER_VIDEO_IS_IN} {VIDEO_FILENAME} {NUMBER OF MINUTES IN VIDEO (integer)}
video_path=$1
video_name=$2;
video_length_minutes=$3
trimmed_name="$(echo $video_name | tr -d '[:space:]')"
trimmed_name="${trimmed_name%.*}"
echo $trimmed_name
mkdir -p images


for ((i=1; i<=$video_length_minutes; i++)) ; do outpath=$video_path/images/${trimmed_name}_$i.jpg && echo $outpath && ffmpeg -ss `echo $i*60.0 | bc` -i "$video_path/$video_name"  -frames:v 1 $outpath; done


