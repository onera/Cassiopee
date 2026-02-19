/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../Data.h"

#ifndef INT64_C
#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)
#endif

#define BYTE uint8_t
#define sR (BYTE)(buffer[s])
#define sG (BYTE)(buffer[s+1])
#define sB (BYTE)(buffer[s+2])

#ifdef HAVE_AV_CONFIG_H
#undef HAVE_AV_CONFIG_H
#endif

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavutil/mem.h>
//#include <libavutil/opt.h>
#include <libavutil/imgutils.h>
}

//=============================================================================
/*
  Write buffer to MPEG file.
  si mode=1, ferme le fichier et nettoie.
*/
//=============================================================================
void writeMPEGFrame(Data* d, char *filename, char* buffer, 
                    E_Int width, E_Int height, E_Int mode)
{
  int gotOutput, ret;
  AVPacket pkt;

  // Initilisation des codecs et du context
  AVCodecContext* context = NULL;
  if (d->ptrState->context == NULL)
  { 
    avcodec_register_all();

    /* find the video encoder */
    AVCodecID codecId1 = AV_CODEC_ID_MPEG1VIDEO; // OK
    AVCodecID codecId2 = AV_CODEC_ID_MPEG2VIDEO; // OK
    AVCodecID codecId3 = AV_CODEC_ID_MPEG4; // OK
    AVCodecID codecId4 = AV_CODEC_ID_WMV1; // NOT OK
    AVCodecID codecId5 = AV_CODEC_ID_H264;
    AVCodecID codecId; 
    AVCodec *codec;
   
    // Essaie de trouver le meilleur codec
    /*
    codecId = codecId5;
    codec = avcodec_find_encoder(codecId);
    if (!codec)
    {
      codecId = codecId3;
      codec = avcodec_find_encoder(codecId);
      if (!codec)
      {
        codecId = codecId2;
        codec = avcodec_find_encoder(codecId);
        if (!codec)
        {
          codecId = codecId1;
          codec = avcodec_find_encoder(codecId);
          if (!codec) 
          {
            fprintf(stderr, "Error: CPlot: no available codec.\n");
            return;      
          }
        }
      }
    }
    */
    // DBX - impose le codec MPEG2
    codecId = codecId2;
    codec = avcodec_find_encoder(codecId);
    // END DBX

    context = avcodec_alloc_context3(codec);
    d->ptrState->context = (void*)context;

    /* put sample parameters */
    context->bit_rate = 10*width*height; // + eleve -> + qualite
    /* resolution must be a multiple of two */
    context->width = width;
    context->height = height;
    /* frames per second */
    context->time_base = (AVRational){1,25};
    context->gop_size = 250; /* emit one intra frame every X frames */
    context->max_b_frames = 1; // must be 1
#if LIBAVCODEC_VERSION_MAJOR > 56
    context->pix_fmt = AV_PIX_FMT_YUV420P;
#else
    context->pix_fmt = PIX_FMT_YUV420P;
#endif
    //context->compression_level = 9;

    // Set profile to baseline
    //if (codecId == AV_CODEC_ID_H264) av_opt_set(context->priv_data, "preset", "slow", 0);
    //if (codecId == AV_CODEC_ID_MPEG4) av_opt_set(context->priv_data, "preset", "slow", 0);

    if (avcodec_open2(context, codec, NULL) < 0) 
    {
      fprintf(stderr, "Error: CPlot: could not open codec.\n");
      return;
    }
    if (codecId == codecId1) 
      printf("Info: CPlot: Init codecs [MPEG1].\n");
    else if (codecId == codecId2) 
      printf("Info: CPlot: Init codecs [MPEG2].\n");
    else if (codecId == codecId3) 
      printf("Info: CPlot: Init codecs [MPEG4].\n");
    else if (codecId == codecId5) 
      printf("Info: CPlot: Init codecs [H.264].\n");
  }
  else
  {
    context = (AVCodecContext*)d->ptrState->context;
  }
    
  /* ouverture du fichier */
  FILE *f;
  if (d->ptrState->ptrFile == NULL)
  {
    f = d->fopenw(filename, "wb");
    if (!f) 
    {
      fprintf(stderr, "Error: CPlot: could not open file %s.\n", filename);
      return;
    }
    d->ptrState->ptrFile = f;
  }
  else f = d->ptrState->ptrFile;
  
  if (mode == 1) // finalize
  {
    uint8_t endcode[] = { 0, 0, 1, 0xb7 };
    av_init_packet(&pkt);
    pkt.data = NULL; // packet data will be allocated by the encoder
    pkt.size = 0;
    gotOutput = 1;
    while (gotOutput != 0)
    {
      fflush(stdout);
      ret = avcodec_encode_video2(context, &pkt, NULL, &gotOutput);
      if (ret < 0) 
      {
        fprintf(stderr, "Error: CPlot: error encoding frame.\n");
        return;
      }
      if (gotOutput) 
      {
        printf("Write frame (size=%5d)\n", pkt.size);
        fwrite(pkt.data, 1, pkt.size, f);
        av_free_packet(&pkt);
      }
    }
    
    /* add sequence end code to have a real mpeg file */
    fwrite(endcode, 1, sizeof(endcode), f);
    fclose(f);
    avcodec_close(context);
    av_free(context);
    printf("File %s closed.\n", filename);
    d->ptrState->shootScreen = 0;
    d->ptrState->continuousExport = 0;
    return; 
  }

  /* frame */

#if LIBAVCODEC_VERSION_MAJOR > 56
  AVFrame* frame = av_frame_alloc();
  frame = av_frame_alloc();

#else
  AVFrame* frame = avcodec_alloc_frame();
  frame = avcodec_alloc_frame();
#endif

  if (!frame) 
  {
    fprintf(stderr, "Error: CPlot: could not allocate video frame.\n");
    return;
  }
  frame->format = context->pix_fmt;
  frame->width = context->width;
  frame->height = context->height;

  ret = av_image_alloc(frame->data, frame->linesize, 
                       context->width, context->height,
                       context->pix_fmt, 32);
  if (ret < 0) 
  {
    fprintf(stderr, "Error: CPlot: could not allocate raw picture buffer.\n");
    return;
  }

//   /* alloc image and output buffer */
//   size = context->width * context->height;
//   uint8_t* frameBuf = (uint8_t*)malloc( ((size*3)/2)*sizeof(uint8_t)); /* size for YUV 420 */
//   frame->data[0] = frameBuf;
//   frame->data[1] = frame->data[0] + size;
//   frame->data[2] = frame->data[1] + size / 4;
//   frame->linesize[0] = context->width;
//   frame->linesize[1] = context->width / 2;
//   frame->linesize[2] = context->width / 2;


  /* passage en YUV420P */
  //unsigned int i = 0;
  unsigned int s = 0;
  int jp, kp, jm;
     
  for (int j = 0; j < height; j++)
    for (int k = 0; k < width; k++)
    {
      jm = height-1-j;
      // Y
      frame->data[0][jm*frame->linesize[0]+k] = (BYTE)( (66*sR + 129*sG + 25*sB + 128) >> 8) + 16;
          
      // UV
      if (0 == j%2 && 0 == k%2)
      {
        jp = jm/2; kp = k/2;
        frame->data[1][jp*frame->linesize[1]+kp] = (BYTE)( (-38*sR - 74*sG + 112*sB + 128) >> 8) + 128;
        frame->data[2][jp*frame->linesize[2]+kp] = (BYTE)( (112*sR - 94*sG - 18*sB + 128) >> 8) + 128;
      }
      s += 3;
    }

  // encode
  av_init_packet(&pkt);
  pkt.data = NULL; // packet data will be allocated by the encoder
  pkt.size = 0;
  frame->pts = d->ptrState->exportNumber;
  ret = avcodec_encode_video2(context, &pkt, frame, &gotOutput);
  if (ret < 0) {
    fprintf(stderr, "Error: CPlot: error encoding frame.\n");
    return;
  }
  if (gotOutput) 
  {
    printf("Write frame (size=%5d)\n", pkt.size);
    fwrite(pkt.data, 1, pkt.size, f);
    av_free_packet(&pkt);
  }
  
  av_freep(&frame->data[0]);
#if LIBAVCODEC_VERSION_MAJOR > 56
  av_frame_free(&frame);
#else
  avcodec_free_frame(&frame);
#endif
}
