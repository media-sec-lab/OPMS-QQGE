% -------------------------------------------------------------------------
% Copyright (c) 2020
% Rémi Cogranne, UTT (Troyes University of Technology)
% All Rights Reserved.
% -------------------------------------------------------------------------
% This code is provided by the author under GNU GENERAL PUBLIC LICENSE GPLv3
% which, as explained on this webpage
% https://www.gnu.org/licenses/quick-guide-gplv3.html
% Allows modification, redistribution, provided that:
% * You share your code under the same license ;
% * That you give credits to the authors ;
% * The code is used only for Education or Research purposes.
% -------------------------------------------------------------------------

% Only needed for fair execution time comparisons
maxNumCompThreads(1);

%First run seems to be always much slower ... for a meaningful execution time comparison we do a first pre-heating / dry-run
coverStruct = jpeg_read( [ '/data/ALASKA_50774_QF75.jpg' ] );
[stegoStruct , pChange , ChangeRate , Deflection ] = JMiPOD_fast(coverStruct, 0.2);

imgList = dir('/data/*.jpg');
for imgIdx = 1 : numel(imgList) ,
    fprintf('\n  ***** Processing image %s ***** \n' , imgList(imgIdx).name );
    Payload = 0.4;
    %read DCT coefficients from JPEG file
    coverStruct = jpeg_read( [ imgList(imgIdx).folder '/' imgList(imgIdx).name ] );
    %get stego DCT coefficients (and Deflection, pChanges and overall ChangeRate)
    tStart = tic;
    [stegoStruct , pChange , ChangeRate , Deflection ] = JMiPOD_fast(coverStruct, Payload);
    tEnd = toc(tStart);
    jpeg_write(stegoStruct , [ '/results/' imgList(imgIdx).name ])
    StegoDCT = coverStruct.coef_arrays{1};
    nbnzAC = ( sum( sum( StegoDCT ~= 0) ) - sum( sum( StegoDCT(1:8:end, 1:8:end) ~= 0 ) ) );
    pChange = pChange/2;
    HNats = -2* pChange(:) .* log( pChange(:) ) - (1- 2* pChange(:) ) .* log( 1- 2* pChange(:) ) ;
    Hbits = -2* pChange(:) .* log2(pChange(:) ) - (1- 2* pChange(:) ) .* log2(1- 2* pChange(:) ) ;
	fprintf("\t\t\t\t\t\t\t\t  Target payload = %5.2f bits\n", Payload*nbnzAC )
	fprintf("JMiPOD fast implemetations runs in %2.3f sec. \t\t\tActual payload  :  %5.2f bits = %5.2f Nats (binary entropy computed from pChanges) \n" , tEnd , nansum( Hbits ) , nansum( HNats ) )

    tStart = tic;
    [stegoStruct , pChange , ChangeRate , Deflection ] = JMiPOD(coverStruct, Payload);
    tEnd = toc(tStart);
    pChange = pChange/2;
	fprintf("JMiPOD original (not so fast) runs in %2.3f sec.\t\tActual payload  :  %5.2f bits = %5.2f Nats (binary entropy computed from pChanges) \n" , tEnd , nansum( -2* pChange(:) .* log2( pChange(:) ) -  (1- 2* pChange(:) )  .* log2(1- 2* pChange(:) ) ) , nansum( -2* pChange(:) .* log( pChange(:) ) -  (1- 2* pChange(:) )  .* log(1- 2* pChange(:) ) ) )

    tStart = tic;
    [ S_STRUCT , pChange , ChangeRate ] = UERD(coverStruct, Payload);
    tEnd = toc(tStart);
    pChange = pChange/2;
	fprintf("UERD runs in %2.3f sec.\t\t\t\t\t\tActual payload  :  %5.2f bits = %5.2f Nats (binary entropy computed from pChanges) \n" , tEnd , nansum( -2* pChange(:) .* log2( pChange(:) ) -  (1- 2* pChange(:) )  .* log2(1- 2* pChange(:) ) ) , nansum( -2* pChange(:) .* log( pChange(:) ) -  (1- 2* pChange(:) )  .* log(1- 2* pChange(:) ) ) )

    coverImage = double( imread( [ imgList(imgIdx).folder '/' imgList(imgIdx).name ] ) );
    tStart = tic;
    [ S_STRUCT , pChange , ChangeRate ] = J_UNIWARD(coverStruct , coverImage , Payload);
    tEnd = toc(tStart);
    pChange = pChange/2;
	fprintf("JUNIWARD runs in %2.3f sec.\t\t\t\t\tActual payload  :  %5.2f bits = %5.2f Nats (binary entropy computed from pChanges) \n" , tEnd , nansum( -2* pChange(:) .* log2( pChange(:) ) -  (1- 2* pChange(:) )  .* log2(1- 2* pChange(:) ) ) , nansum( -2* pChange(:) .* log( pChange(:) ) -  (1- 2* pChange(:) )  .* log(1- 2* pChange(:) ) ) )
end