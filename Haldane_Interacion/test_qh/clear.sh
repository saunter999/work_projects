if [ -n "$1" ]
then
  if [ $1 = "-d" ] || [ $1 = "--diag" ]; then
      rm -rf ./diagram/*
  elif [ $1 = "-a" ] || [ $1 = "--all" ]; then
      rm *.hkl
      rm *.txt
      rm statis_total.hkl
      rm _job*.sh
  fi
fi
rm dyson/Message.txt
rm dyson/*.hkl
rm dyson/*.txt
rm dyson/*.pdf
#rm Coordinates.txt
rm *.log
rm *.gv
rm -rf infile
rm -rf outfile
#rm -rf ./diagram/
#mkdir ./diagram
