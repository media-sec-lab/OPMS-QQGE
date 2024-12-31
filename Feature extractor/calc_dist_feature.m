function dist=calc_dist_feature(cover_path,stego_path,quality,fea_name)

if strcmp(fea_name,'DCTR')==1
    fea_c = DCTR({cover_path},quality);
    fea_s = DCTR({stego_path},quality);
elseif strcmp(fea_name,'GFR')==1
    fea_c = GFR(cover_path,32,quality);
    fea_s = GFR(stego_path,32,quality);
end
dist=sqrt(sum((fea_c-fea_s).^2));

end