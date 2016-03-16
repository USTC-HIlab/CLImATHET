function scores = Scale_score(scores,m_score)

wf = 0.7;
score_th = wf*m_score;
scores = (scores-score_th)/score_th+scores;
scores(scores > 1) = 1;
scores(scores < 0) = 0;
scores = scores*100;

end