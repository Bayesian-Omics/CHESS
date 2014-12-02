/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from xml_file_write.cc in the ESS++ program
 *      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
 *                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
 *                    David Hastie (d.hastie@imperial.ac.uk)
 *      2010
 *
 * The file is modified from xml_file_write.cc in the FREGENE program
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
#include <string.h>
#include <time.h>
#include <math.h>
#include "xml_file_write.h"

MaXmlTagWrite::MaXmlTagWrite()
{
	XF=NULL;
}

MaXmlTagWrite::MaXmlTagWrite(FILE *f, char version[10])
{
	SetXmlFile(f);
	WriteDeclaration(version);
}

MaXmlTagWrite::MaXmlTagWrite(char *s, char version[10])
{
	XF=OpenXmlFile(s);
	WriteDeclaration(version);
}

MaXmlTagWrite::~MaXmlTagWrite()
{
	if(XF) fclose(XF);
	XF=NULL;
}

FILE* MaXmlTagWrite::OpenXmlFile(char *s)
{
	FILE *fx=NULL;

	fx=fopen(s,"wb");

	return fx;
}

int MaXmlTagWrite::CloseXmlFile(FILE *f)
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

void MaXmlTagWrite::SetXmlFile(FILE *f) 
{
	if(f) XF=f;
	else{
	  printf("usage:: common/xml_file_read.cc\n");
	  printf("Can not set XML file\n");
	  exit(1);
	}
}

int MaXmlTagWrite::WriteDeclaration(char version[10])
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"<?xml version='%s' ?>\n",version);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::WriteDocType()  //to perform
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"\n");
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::WriteComment(char *s)  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"<!-- %s -->",s);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;

	return n;
}

int MaXmlTagWrite::WriteString(char *s)  
{	
	if(!XF) return MA_XML_FILE_ERROR;
	int tn =strlen(s);
	char *command=new char[tn+1];
	sprintf(command,"%s",s);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;

	return n;
}

int MaXmlTagWrite::WriteSpace()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	int n=fwrite("\n",1,strlen("\n"),XF);

	return n;
}

int MaXmlTagWrite::WriteLine(char *s)  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"%s",s);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;

	return n;
}

int MaXmlTagWrite::WriteTab()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	int n=fwrite("\t",1,strlen("\t"),XF);

	return n;
}

int MaXmlTagWrite::WriteEmptyElement(char *s)
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"<%s/>\n",s);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::StartTag(char *s)
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"<%s>",s);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::EndTag(char *s)
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"</%s>\n",s);
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::StartCDATA()
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"<![CDATA[\n");
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::EndCDATA()
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	sprintf(command,"\n]]>\n");
	int n=fwrite(command,1,strlen(command),XF);

	delete [] command;
	return n;
}

int MaXmlTagWrite::StartHead()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	StartTag("HEAD");

	return MA_NO_ERROR;
}

int MaXmlTagWrite::EndHead()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	EndTag("HEAD");

	return MA_NO_ERROR;
}

int MaXmlTagWrite::StartBinData()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	StartTag("DATA");
	StartCDATA();

	return MA_NO_ERROR;
}

int MaXmlTagWrite::EndBinData()  
{	
	if(!XF) return MA_XML_FILE_ERROR;
	
	EndCDATA();
	EndTag("DATA");

	return MA_NO_ERROR;
}

int MaXmlTagWrite::StartAData()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	StartTag("DATA");

	return MA_NO_ERROR;
}

int MaXmlTagWrite::EndAData()  
{	
	if(!XF) return MA_XML_FILE_ERROR;
	
	EndTag("DATA");

	return MA_NO_ERROR;
}

int MaXmlTagWrite::WriteAuthor(char *s)  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	StartTag("AUTHOR");
	int n=fwrite(s,1,strlen(s),XF);
	EndTag("AUTHOR");

	return n;
}

int MaXmlTagWrite::WriteDate()  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	StartTag("DATE");
   time_t ltime;
   time( &ltime );
   strcpy(command,ctime( &ltime ));
	int n=fwrite(command,1,strlen(command),XF);
	EndTag("DATE");

	delete [] command;
	return n;
}

int MaXmlTagWrite::WriteType(char *s)  
{	
	if(!XF) return MA_XML_FILE_ERROR;

	char *command=new char[1024];
	StartTag("TYPE");
	if((strcmp("T",s)!=0) || (strcmp("B",s)!=0)) strcpy(command,"UT");
	else strcpy(command,s);
	int n=fwrite(command,1,strlen(command),XF);
	EndTag("TYPE");

	delete [] command;
	return n;
}



