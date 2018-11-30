filtered_out_cats <- "A" 

get("GO:0051171", GOTERM)
Ontology("GO:0051171")
ancestors <- get("GO:0051171", GOBPANCESTOR)
for (ancestor in ancestors) {
  print("==========")
  print(get(ancestor, GOTERM))
}

term <- "GO:0005634"
get(term, GOTERM)
Ontology(term)

ancestors <- get(term, GOCCANCESTOR)
for (ancestor in ancestors) {
  print("==========")
  print(get(ancestor, GOTERM))
}


term <- "GO:0005634"
offsprings <- get(term, GOCCOFFSPRING)
print(offsprings)
for (offspring in offsprings) {
  print("==========")
  print(get(offspring, GOTERM))
}