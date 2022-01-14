#!bash
for i in *.tif; do sips -s format png $i --out ${i%.*}.png; done
