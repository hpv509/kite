#!/bin/bash
rsync -avhzp --exclude 'build/' ~/devel/fork/sync/zv/ "viking:~/scratch/Tmp/"
