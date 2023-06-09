# wetland_bird_acoustic

Welcome to the project repository for my master's thesis: What can the
soundscape of ephemeral wetlands tell us about the drivers of songbird community
composition and diversity?

## Objective
The objective of this project is to collect acoustic samples of songbirds around
ephemeral wetlands in the Lowcountry of South Carolina and compare community
composition results to the results produced by traditional point count methods.

## Background/Question/Methods
Wetlands provide essential ecosystem services and natural habitats for many
organisms, however, ephemeral wetlands are often overlooked because it is
unclear if these partially wet areas are ecologically unique and important.
Traditional survey methods provide essential data that informs the diversity of
many songbird communities, however, they are often costly and cumbersome when
compared to passive acoustic monitoring methods. Our study characterizes the
soundscape of ephemeral wetlands to examine whether automated methods are
redundant or complementary to human detection methods. We compare community
diversity and composition of songbirds in two ways: 1) in wet vs. dry sites and
2) using human vs. automated detection. The soundscape is characterized using
solo system audio recorders and the BirdNET artificial intelligence algorithm to
identify species occurrences. We sampled 24 isolated forest ephemeral wetlands
and 6 open-canopy longleaf pine savanna uplands for 3 days each between May 15 -
June 15, 2022. We compare the automated species detections to traditional
point-count surveys to inform future monitoring schemes.

## Project Structure

The structure of this code-base is R. To access and understand the relevant data,
download the following files:

./data/compiled_bird_audio_master.csv
./data/metadata.csv

To recreate analysis run the following scripts in order:

./scripts/data_processing.R
./scripts/analysis.R

My results can be replicated by using the various plot functions
included in the analysis code to illustrate the relationship between different
variables. I utilize linear models with assumed normality.

## Acknowledgements

I would like to thank Dr. Dan McGlinn for being the best advisor I could ask for
and for his inspiring curiosity towards ecological questions and processes. I
appreciate Jackson Barratt Heitmann's leadership in the field and Sean Cannon's
assistance collecting vegetation data.
