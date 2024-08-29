graph [
  name "GO tree"
  graph [
  ]
  node [
    id 0
    label "GO:0110165\ncellular anatomical entity"
    fillcolor "plum"
  ]
  node [
    id 1
    label "GO:0005575\ncellular_component"
    fillcolor "plum"
  ]
  node [
    id 2
    label "GO:0000323\nlytic vacuole"
    fillcolor "plum"
  ]
  node [
    id 3
    label "GO:0005764\nlysosome"
    fillcolor "plum"
  ]
  node [
    id 4
    label "GO:0005773\nvacuole"
    fillcolor "plum"
  ]
  node [
    id 5
    label "GO:0043229\nintracellular organelle"
    fillcolor "plum"
  ]
  node [
    id 6
    label "GO:0043226\norganelle"
    fillcolor "plum"
  ]
  node [
    id 7
    label "GO:0043231\nintracellular membrane-bounded organelle"
    fillcolor "plum"
  ]
  node [
    id 8
    label "GO:0043227\nmembrane-bounded organelle"
  ]
  edge [
    source 0
    target 6
  ]
  edge [
    source 0
    target 1
  ]
  edge [
    source 2
    target 3
  ]
  edge [
    source 2
    target 4
  ]
  edge [
    source 4
    target 7
  ]
  edge [
    source 5
    target 7
  ]
  edge [
    source 5
    target 6
  ]
  edge [
    source 6
    target 8
  ]
  edge [
    source 7
    target 8
  ]
]
