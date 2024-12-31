1. Please refer to "OPMS_QQGE_main.m" for the main code.

2. Notes on using some deblocking models:
- SSRQC: Official code, located in the "SSRQC" folder.  
- DnCNN: Models retrained using the official MATLAB version code, located in the "DnCNN" folder.  
- DRUNet and SwinIR: Models fine-tuned on the official pre-trained parameters, located in the "drunet.zip" and "swinir.zip" archives. Images of different QFs utilize their respective trained deblocking models. 
   * DRUNet
      Model for images of QF50和QF75: ...\drunet\DPIR-master\model_zoo\drunet_deblocking_grayscale.pth
      Model for images of QF95: ...\drunet\DPIR-master\model_zoo\50000_G.pth
   * SwinIR
      Model for images of QF50: ...\SwinIR\model_zoo\swinir\dejpeg_50\006_CAR_DFWB_s126w7_SwinIR-M_jpeg40.pth
      Model for images of QF75: ...\SwinIR\model_zoo\swinir\dejpeg_75\80000_G.pth
      Model for images of QF95: ...\SwinIR\model_zoo\swinir\dejpeg_95\30000_G.pth
