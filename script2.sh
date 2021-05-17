counter2=69
counter=1
counter3=16
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib pruebawt.cpp -o program4 -lsdsl
while [ $counter -le 5 ]
do
string1="salida"$counter
g++ generador.cpp
./a.out  tiger_map_usa2.hierarchy > salidaprueba$counter
./program4 tiger_map_usa_filtered_6lvls.pg tiger_map_usa_filtered_6lvls.hierarchy salidaprueba$counter nuevosresultados$counter2 tamP$counter2 vecinosP$counter2 $counter3  > salidapropio$counter
echo $string1
((counter3=counter3*2))
((counter++))
done
