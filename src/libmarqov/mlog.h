/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2022, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MLOG_H
#define MLOG_H

#define FLOG_DISABLE_VERBOSE_FUN
#include "flog/flog.h"

#define MLOG_EXTLEVEL FLOG_EXTLEVEL

#define MLOGDEBUGVERBOSE FLOGDEBUGVERBOSE
#define MLOGDEBUG FLOGDEBUG
#define MLOGRELEASEVERBOSE FLOGRELEASEVERBOSE
#define MLOGVERBOSE FLOGVERBOSE
#define MLOGRELEASE FLOGRELEASE

#endif
