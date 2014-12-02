/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from xml_file_read.cc in the ESS++ program
 *      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
 *                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
 *                    David Hastie (d.hastie@imperial.ac.uk)
 *      2010
 *
 * The file is modified from xml_file_read.cc in the FREGENE program
 *      Copyright (c) Clive J Hoggart   <c.hoggart@imperial.ac.uk>
 *                    Taane G Clark     <tgc@well.ox.ac.uk>
 *                    Marc Chadeau      <m.chadeau@imperial.ac.uk>
 *
 * CHESS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CHESS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CHESS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "xml_file_read.h"

MaXmlTagRead::MaXmlTagRead(FILE *f)
{
	SetXmlFile(f);
}

MaXmlTagRead::MaXmlTagRead(char *s)
{
	XF=OpenXmlFile(s);
}

MaXmlTagRead::~MaXmlTagRead()
{
	//if(XF) fclose(XF);
	//XF=NULL;
}

void MaXmlTagRead::SetXmlFile(FILE *f) 
{
	if(f) XF=f;
	else{
	  printf("usage:: common/xml_file_read.cc\n");
	  printf("Can not set XML file\n");
	  exit(1);
	}
}

FILE* MaXmlTagRead::OpenXmlFile(char *s)
{
	FILE *fx=NULL;

	fx=fopen(s,"rb");

	return fx;
}

int MaXmlTagRead::CloseXmlFile(FILE *f)
{
	if(f) 
		if(fclose(f)) return MA_CLOSE_FILE_ERROR;
		else 
		{
			f=NULL;
			return MA_NO_ERROR;
		}
   else return MA_NO_ACTION;
}

int MaXmlTagRead::CloseXmlFile()
{
	if(XF) 
		if(fclose(XF)) return MA_CLOSE_FILE_ERROR;
		else 
		{
			XF=NULL;
			return MA_NO_ERROR;
		}
   else return MA_NO_ACTION;
}

//*********************************************************
//*	PURPOSE: 
//*		
//*	IN: tag name
//*
//*	OUTPUT: 
//*
//*	RETURN: Return true if found else false
//********************************************************
bool MaXmlTagRead::GetTagPos(const char *tag, long start, long end, long &pos_found)
{
   long notused=end;
   end=notused;
	bool found = false;
	bool endtag= false;

	if(strcmp(tag,">")==0) endtag=true;
	char *tab = new char[256];
	memset(tab,0,256*sizeof(char));
	
	int nn=strlen(tag);

	int ii=1, NN;

	if(!endtag) 
	{
		tab[0]='<';
		for(int i=0; i<nn;i++) tab[i+1]=tag[i];
		tab[nn+1]='\0';
		NN=nn+1;
	}
	else 
	{
		for(int i=0; i<nn;i++) tab[i]=tag[i];
		tab[nn]='\0';
		NN=nn;
	}
	
	fseek(XF,start,SEEK_SET);
	long CURRENT_POS=start;
	char ch;

	while(!found && !feof(XF))
	{
		ii=0;

		while(!feof(XF) && (fread(&ch,sizeof(char),1,XF) == 1))
		{
			if(tab[ii]==ch) 
			{
				ii++;
			}
			else ii=0;

			if(ii==NN)
			{
				CURRENT_POS=ftell(XF);
				pos_found=CURRENT_POS;
				found = true;
				break;
			}
		}

		if(found && !endtag) 
		{
			long pos;
			pos=ftell(XF);
			fread(&ch,sizeof(char),1,XF);
			fseek(XF,pos,SEEK_SET);

			if(ch=='>' || ch==' ') GetTagPos(">", CURRENT_POS, 0, pos_found);
			else found=false;
		}
	}
	

	delete [] tab;
	return found;
}

//********************************************************************
//*	PURPOSE: 
//*		
//*	IN: tag name
//*
//*	OUTPUT: the first position before the last end tag "tag name"
//*
//*	RETURN: Return true if found else false
//********************************************************************
bool MaXmlTagRead::GetTagFirstPos(const char *tag, long start, long end, long &first_pos)
{
   long notused=end;
   end=notused;
	bool found = false;
	bool endtag= false;

	if(strcmp(tag,">")==0) endtag=true;
	char *tab = new char[256];
	memset(tab,0,256*sizeof(char));

	int nn=strlen(tag);

	int ii=1, NN;
	char CH1=tag[0];
	
	if(!endtag) 
	{
		tab[0]='<';
		for(int i=0; i<nn;i++) tab[i+1]=tag[i];
		tab[nn+1]='\0';
		NN=nn+1;
	}
	else 
	{
		for(int i=0; i<nn;i++) tab[i]=tag[i];
		tab[nn]='\0';
		NN=nn;
	}

	
	fseek(XF,start,SEEK_SET);
	long CURRENT_POS=start;
	char ch;


	while(!found && !feof(XF))
	{
	
		ii=0;

		while(!feof(XF) && (fread(&ch,sizeof(char),1,XF) == 1))
		{
			if(tab[ii]==ch) 
			{
				ii++;
			}
			else ii=0;

			if(ch==CH1) 
			{
				first_pos=ftell(XF);
				first_pos-=3;
			}
			if(ii==NN)
			{
				CURRENT_POS=ftell(XF);
				found = true;
				break;
			}

		}

		if(found && !endtag) 
		{
			long pos;
			pos=ftell(XF);
			fread(&ch,sizeof(char),1,XF);
			fseek(XF,pos,SEEK_SET);
			long pos_found;
			if(ch=='>' || ch==' ') GetTagPos(">", CURRENT_POS, 0, pos_found); //pourra servir pour tester l'excistance de >
			else found=false;
		}
	}

	delete [] tab;

	return found;
}

//*********************************************************
//*	PURPOSE: read the tag content
//*		
//*	IN: tag name
//*
//*	OUTPUT: the content between tag start and tag end
//*
//*	RETURN: Return true if found else false
//********************************************************
bool MaXmlTagRead::ReadTag(const char *tag, long start, long end, char *content, unsigned long content_size)
{
	long pos_found, first_end;

	if(!GetTagPos(tag, start, end, pos_found)) return false;

	char *endtag = new char[256];
	sprintf(endtag,"/%s",tag);
	if(!GetTagFirstPos(endtag, start, end, first_end))
	{
		delete [] endtag;
		return false;
	}

	fseek(XF,pos_found,SEEK_SET);
	unsigned long N=(unsigned int)(first_end-pos_found+1);

	N=Mymin(content_size,N);

	if(fread(content,sizeof(char),N,XF) != N) 
	{
		delete [] endtag;
		return false;
	}

	content[N]='\0';

	delete [] endtag;
	return true;
}

bool MaXmlTagRead::ReadTag(const char *tag, long start, long end, char *content, unsigned long content_size, long &pos_found)
{
	long first_end;

	if(!GetTagPos(tag, start, end, pos_found)) return false;

	char *endtag = new char[256];
	sprintf(endtag,"/%s",tag);
	if(!GetTagFirstPos(endtag, start, end, first_end))
	{
		delete [] endtag;
		return false;
	}

	fseek(XF,pos_found,SEEK_SET);
	unsigned long N=(unsigned int)(first_end-pos_found+1);

	N=Mymin(content_size,N);

	if(fread(content,sizeof(char),N,XF) != N) 
	{
		delete [] endtag;
		return false;
	}

	content[N]='\0';

	delete [] endtag;
	return true;
}

//*********************************************************
//*	PURPOSE: read the size of the tag content
//*		
//*	IN: tag name
//*
//*	OUTPUT: the size of the tag content
//*
//*	RETURN: Return true if found else false
//********************************************************
bool MaXmlTagRead::ReadSizeTag(char *tag, long start, long end, unsigned long &content_size)
{
	long pos_found, first_end;

	if(!GetTagPos(tag, start, end, pos_found)) return false;

	char *endtag = new char[256];
	sprintf(endtag,"/%s",tag);
	if(!GetTagFirstPos(endtag, start, end, first_end))
	{
		delete [] endtag;
		return false;
	}

	fseek(XF,pos_found,SEEK_SET);
	unsigned long N=(unsigned int)(first_end-pos_found+1);

	content_size=N;

	delete [] endtag;
	return true;
}

//*********************************************************
//*	PURPOSE: 
//*		
//*	IN: 
//*
//*	OUTPUT: 
//*
//*	RETURN: 
//********************************************************
char *MaXmlTagRead::GetTagEndOf(char *tag)
{

	char *endtag = new char[256];

	sprintf(endtag,"/%s",tag);

	return endtag;

}
