# LOOPS as in python

animals = ["dog", "cat", "mouse"]

for animal in animals
    println("$animal is a mammal")
end

# with enumeration
for (i,animal) in enumerate(animals)
    println("$animal no. $i is a mammal")
end

verbs = ["plays","sleeps","eats"]
for (animal,does) in zip(animals,verbs)
    println("The $animal $does.")
end
