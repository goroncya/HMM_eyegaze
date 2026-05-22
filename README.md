# HMM_eyegaze
HMM study on the impact of skin abberations on facial exploration behavior

The R code used to fit the Participant x Media HMMs, perform feature extraction, implement state matching, and conduct the clustering analyses consists of the following R scripts:
1.	build_HMM.R, defining functions:
  a.  build_hmm_multiseq()
  b.	get_state_table()
  c.	get_state_mapping()
  d.	extract_state_params()
  e.	extract_transition_vector()
  f.	get_state_info()
  g.	calc_state_pair_metrics()
2. Media_comparison_1.R
3. Media_comparison_2.R
4. Media_comparison_3.R
5. Media_comparison_4.R

Input data structure
The analysis is based on a data frame named 'df', where each row corresponds to a single gaze observation (fixation). The data frame contains the following columns:
- 'ParticipantID' - unique identifier of the participant,
- 'MediaID' - identifier of the viewed stimulus/image,
- 'ModelID' - identifier of the person/model presented within a given media item. In the HMM analysis, each 'ModelID' is treated as a separate observation sequence within a given participant–media pair,
- 'TrialID' - order of the observation within the sequence,
- 'OX' - horizontal gaze coordinate,
- 'OY' - vertical gaze coordinate.
  
References: Marek Jankowski, Krzysztof Jasiński, Agnieszka Goroncy, Impact of skin abberations on facial exploration behavior: an HMM study, under review, 2026.

