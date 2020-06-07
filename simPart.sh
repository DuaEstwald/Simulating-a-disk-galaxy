echo '#\t Npart time\n' > timesimulation.txt
echo '#Without Self-gravity\n' >> timesimulation.txt

for n in 10 100 1000
do
	python Part2.py $n
done

echo '#Self-gravity\n' >> timesimulation.txt

for n in 10 100 1000
do
	python Part3.py $n
done

echo '#Self-gravity + 10KDTree\n' >> timesimulation.txt

for n in 10 100 1000
do 
	python EXTRApart.py $n 10
done

echo '#Self-gravity + 50KDTree\n' >> timesimulation.txt

for n in 100 1000
do 
	python EXTRApart.py $n 50
done


