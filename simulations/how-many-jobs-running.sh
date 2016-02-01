#!/bin/bash

#COMPUTERS="garnish juno"
#COMPUTERS="garnish"
#COMPUTERS="equinox"

COMPUTERS="buffoon cave equinox garnish glencoe goldrush grumpy haystack hiccup nordic norinori nosey otter oval pentopia platypus primrose rainbow saratoga senate sloth snake solstice stadium suraromu tomtom tapa"

FASTCOMPUTERS="gargoyle orwell snake garnish chiffchaff heraclius jigsaw norinori eris numberlink pluto janus minesweeper chopin haydn shikaku himalia metis anice statuepark rossini mozart socrates elara thebe sinope brahms paaliaq persaeus mustard magnets wagner scramble beans wryneck basil sun jupiter rubik democrates suguru hashi porridge tangram siarnaq verdi epicurus cayenne"

ALLCOMPUTERS="equinox parrot bonxie buzzard blackcap parakeet gull kite plover whinchat skua sanderling primrose saratoga glencoe treecreeper minesweeper porridge hashi hippo democrates chopin yogi haydn vulture tacitus tapa platypus himalia metis magnets verdi grumpy jigsaw siarnaq nosey norinori hiccup otter eris heraclius numberlink thales pluto tangram janus wagner buffoon acrion mustard thebe persaeus paaliaq neverland sinope brahms rubik anice sloth statuepark rossini tapir jupiter sun basil mozart socrates fulmar epicurus elara inglorious beans scramble diogenes cayenne garnish gargoyle suguru shikaku wryneck stonechat magpie waxwing bunting snake chiffchaff senate cave nordic stadium pentopia haystack oval tomtom quail harrier orwell diver"
#COMPUTERS="buffoon cave demonfan equinox garnish glencoe goldrush grumpy hashi haystack hiccup jigsaw jubilee neverland nordic norinori nosey otter oval pentopia platypus primrose rainbow rubik saratoga senate shikaku sloth snake solstice stadium suraromu tangram tapa tomtom yogi zodiac"
USR="st624"




#  ### wake computers
#  for computer in $COMPUTERS
#  do
#   	(ssh -qf juno wake ${computer})
#  done
#  echo "Waiting for computers to wake up"
#  sleep 2m
#  echo "Finished waiting"


# ## count number of cores
# totproc=0
# for computer in $COMPUTERS
# do
#     echo "$computer"
#     nproc=$(ssh -f $computer "nproc")
#     if [ -n "$nproc" ]; then
#         totproc=$(($totproc+$nproc))
#     fi
# done
# echo "$totproc processors"




### connect to each computer
#/bin/nice -n 19

##### print how many processors I am not using on other computers
#for computer in $COMPUTERS
#do
#    nproc=$(ssh -f $computer "nproc")
#    nrunning=$(ssh $computer "ps -u st624 | grep R | wc -l")
#    if [ -n "$nproc" ]; then  ## only execute the following if the computer responded
#        if [ $(($nproc-1)) -ge $(($nrunning+1)) ]; then
#            echo "$computer"
#            echo "$(($nproc-1-$nrunning))"
#        fi
#    fi
#done



totrunning=0
totproc=0
#wakeup $FASTCOMPUTERS

for computer in $FASTCOMPUTERS
do
    nproc=$(ssh -f $computer "nproc")
    if [ -n "$nproc" ]; then  ## only execute the following if the computer responded 
        nrunning=$(ssh -f $computer "ps -u st624 | grep R | wc -l")
        #        for i in $(seq 1 $(($nproc-1 - $nrunning))); do ## leave one CPU alone
        echo $computer $nrunning
        totrunning=$(($totrunning+$nrunning))
        totproc=$(($totproc+$nproc-1))
    fi
done
echo $totrunning 'jobs running in total'
echo $totproc 'cpus available in total'


############################## remove .Rout if .RData exists

files=$(SIM2016a*.RData)
for file in SIM2016a*.RData
do
 echo "$file"
 echo "${file%.*}"
 rm -v "${file%.*}".Rout
 done

