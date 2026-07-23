# Registro de uso de IA e de prompts

O Prof. Jones Albuquerque solicita que o uso de IA e os prompts empregados sejam
divulgados. Este arquivo descreve, de forma transparente, como a ferramenta de IA
generativa (Anthropic Claude) foi utilizada neste trabalho e qual foi a
participação da autora em cada etapa.

## Como o trabalho foi conduzido

A IA foi usada como ferramenta de apoio dentro desse processo, com pesos
diferentes em cada frente:

- **Texto do artigo:** a IA produziu rascunhos a partir das direções da autora. O
  autor revisou o manuscrito por inteiro, ajustando redação, ênfase e estrutura
  até a versão final, pela qual é responsável.
- **Código:** a implementação foi majoritariamente gerada pela IA a partir das
  especificações da autora (modelo SEIR por Runge-Kutta, autômato celular na
  vizinhança de Moore, contagem de caixas e validação na curva de Koch). Contando com
  ajustes da autora, que executou e conferiu os resultados.
- **Figuras e números:** gerados pela execução do código (`run_experiments.py`) e
  conferidos pela autora antes de entrarem no texto.

## Prompts e direcionamentos principais

Os direcionamentos dados à IA, na ordem em que moldaram o trabalho, foram, em
resumo:

1. Enquadrar a atividade como entrega única para as duas disciplinas, a partir das
   referências indicadas pelo professor (Keeling & Rohani, Cap. 1-2; Schiff,
   Cap. 1; Fuentes & Kuperman, 1999), e aprofundar nelas.
2. Modelar o sarampo com SEIR e desenhar o projeto de modo que a dimensão de
   similaridade fosse a ponte entre o modelo compartimental (EDO) e o autômato
   celular.
3. Implementar o SEIR de campo médio com Runge-Kutta de quarta ordem; implementar
   o autômato celular SEIR 2D na vizinhança de Moore, com núcleo infeccioso
   central; implementar a contagem de caixas e a validação na curva de Koch.
4. Explorar a transmissibilidade para localizar o limiar de invasão e gerar as
   quatro figuras do artigo, salvando os números em `results.json`.
5. Redigir um rascunho do manuscrito em estilo Scientific Reports, em português,
   com base nos números reais gerados pelo código.
6. Ajustes de redação e organização definidos pela autora: título mais formal,
   remoção de travessões e de construções de efeito, inclusão do Abstract em
   inglês, e reordenação para Introdução, Métodos, Resultados e Discussão, com os
   Métodos em texto corrido.

## Responsabilidade

A autora revisou o texto integralmente, executou o código, conferiu os resultados e
assume responsabilidade pelo conteúdo do artigo e do repositório. Nenhuma
referência ou resultado foi inventado: as referências centrais foram indicadas
pelo professor e todos os números remontam a `results.json`.