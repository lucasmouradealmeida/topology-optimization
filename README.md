# **Otimização Topológica**

<p align="justify"> O presente programa foi desenvolvido como parte de uma Iniciação Científica na Universidade Federal do ABC, com orientação do Prof. Dr. Marcelo Araujo da Silva. </p>


## **Introdução**

<p align="justify"> Otimização pode ser definida como um procedimento pelo qual é possível encontrar uma solução ou um conjunto de soluções ótimas para uma determinada função ou um conjunto de funções, que regem um problema específico, submetidas a restrições. Os algoritmos de otimização introduzem estratégias numéricas na busca por soluções ótimas de engenharia. </p>


## **Fundamentação** **Teórica**

<p align="justify"> O procedimento utilizado no presente projeto conta com o uso do Método dos Elementos Finitos (MEF) juntamente com um algoritmo de otimização, para se obter um domínio otimizado. O MEF consiste em discretizar o domínio da estrutura em vários subdomínios, denominados elementos. É montada uma equação para cada elemento e então estas equações são combinadas para se determinar uma expressão que representa a estrutura como um todo, sendo então resolvida e  determinados os deslocamentos, deformações e tensões no domínio da estrutura. Já o algoritmo de otimização consiste em determinar a espessura de cada elemento, de forma a minimizar a função massa total da estrutura, respeitando-se os limites de tensões admissíveis. </p>

<p align="center">
  <img width="401" height="207" src="https://i.imgur.com/M6LU8xf.png">
</p>

<p align="justify"> Na aplicação prática da otimização topológica alguns aspectos são fundamentais. Por exemplo, na Figura acima (topologia obtida) podem surgir pontos com cores intermediárias entre o petro e o branco, denominados de escalas de cinza (ou “gray scale” em inglês). Esses pontos indicam a presença de elementos com espessura intermediária entre a máxima e a nula. Estas espessuras podem não ser viáveis de serem implementadas na prática, mas normalmente ocorrem, ou seja, a presença do “gray scale” é inerente a obtenção da solução ótima. A imagem da estrutura obtida por otimização topológica (OT) representa um excelente  ponto de partida que necessita ser interpretado para se obter o projeto final da estrutura a ser adotado na prática na indústria. Este processo de interpretação se denomina refinamento ou suavização, e pode ser feito utilizando-se métodos de processamento de imagem, ou simplesmente desenhando-se uma estrutura baseada na imagem obtida por OT. Muitas vezes os resultados gerados pela OT não são intuitivos e é necessário fazer uma verificação da estrutura final usando o Método dos Elementos Finitos. A última etapa é a fabricação da estrutura. </p>




