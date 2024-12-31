function  UNIWARD_Sim_Embed(cover_img,stego_img,payload,qt )
%Ä£ÄâÇ¶Èë
[rhoP1,rhoM1]=J_UNIWARD_COST(cover_img,qt);
embedding_3_simulator(cover_img,stego_img,payload,rhoP1,rhoM1);
% embedding_3_stc(cover_img,stego_img,payload,rhoP1,rhoM1);


end

