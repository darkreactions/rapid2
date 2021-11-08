# usage . compress_video.sh {FILENAME OF VIDEO} {PATH_OF_FOLDER_TO_STORE_TEMP_FILES} {HOW MANY FRAMES TO SKIP BETWEEN KEPT FRAMES}
# e.g.:
# . compress_video.sh "MA 333.1a.avi" ~/Downloads/target 1000

export target=$2 ;
export frame_sample_per=$3;
 find . -type f -name "*$1*" \
              -execdir bash -c 'ffmpeg -i "${0}" -vf "select=not(mod(n\,${frame_sample_per}))" \
              -vsync vfr -q:v 2 "${target}/${0}_%03d.png"' {} \; \
              -execdir bash -c 'ffmpeg -r 60 -f image2 -i "${target}/${0}_%03d.png" \
              -vcodec libx264 -crf 25 -pix_fmt yuv420p "${target}/clip_${0}"' {} \; \
              ; rm -f ${target}/*$1_*.png ;
