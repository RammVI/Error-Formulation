/* $Id: blend.h,v 1.2 2000/10/27 15:21:40 mholst Exp $ */

/*
 * Mesa 3-D graphics library
 * Version:  2.0
 * Copyright (C) 1995-1996  Brian Paul
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef BLEND_H
#define BLEND_H


#include "types.h"



extern void gl_blend_span( GLcontext* ctx, GLuint n, GLint x, GLint y,
			   GLubyte red[], GLubyte green[],
			   GLubyte blue[], GLubyte alpha[],
			   GLubyte mask[] );


extern void gl_blend_pixels( GLcontext* ctx,
                             GLuint n, const GLint x[], const GLint y[],
			     GLubyte red[], GLubyte green[],
			     GLubyte blue[], GLubyte alpha[],
			     GLubyte mask[] );

extern void gl_BlendFunc( GLcontext* ctx, GLenum sfactor, GLenum dfactor );

extern void gl_BlendEquation( GLcontext* ctx, GLenum mode );

extern void gl_BlendColor( GLcontext* ctx, GLclampf red, GLclampf green,
                              GLclampf blue, GLclampf alpha );

#endif
