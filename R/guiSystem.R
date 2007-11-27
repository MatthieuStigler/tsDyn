## Copyright (C) 2006  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

#Custom environment S3 class. With methods for subclassing ('extend') and printing
print.MyEnv <- function(x, ...) {
	cat("<",class(x)[1],">\n",sep="")
}

Node <- function(parent=NULL, childs=NULL, opts=NULL, tkvar=NULL, ...) {
	this <- extend(MyEnv(), "Node", parent=NULL, childs=NULL, tkvar=tkvar, opts=opts)
	if(!is.null(parent))
		add(parent, this)
	if(!is.null(childs))
		for(nd in childs)
			add(this, nd)
	return(this)
}

#Add a child to current node
add.Node <- function(this, ...) {
	nodes <- list(...)
	if(length(nodes)>1) {
		for(node in nodes)
			add(this, node)
		return(this)
	}
	this$childs <- c(this$childs, nodes)
	nodes[[1]]$parent <- this
}

#cmd: geometry manager command (as char. string)
#opts: geometry manager options (as char. string)
Frame <- function(cmd="pack", conf=NULL, ...) {
	extend(Node(...), "Frame", cmd=cmd, conf=conf)
}

LabelFrame <- function(cmd="pack", conf=NULL, text="", ...) {
	extend(Frame(cmd="pack", conf=NULL, ...), "LabelFrame", text=text)
}

#opts: widget type and options
Widget <- function(...) {
	extend(Node(...), "Widget")
}

#Creates all necessary (and properly linked) tcltk widgets, and shows them
#'Frame' method
build.Frame <- function(this, ...) {
	parent <- this$parent
	parent.frm <- parent$tkvar
	if(is.null(parent.frm))
		this$tkvar <- tktoplevel()
	else {
		args <- c(list(parent.frm), this$conf)
		this$tkvar <- do.call(tkframe, args)
	}
  nds <- this$childs
  if(!is.null(nds))
  	for(nd in nds)
    	build(nd)
	if(!is.null(parent.frm)) {
  	cmd <- list(parent$cmd, this$tkvar)
  	opts <- parent$opts
  	do.call(tcl, c(cmd, opts))
	}
}

build.LabelFrame <- function(this, ...) {
	parent <- this$parent
	parent.frm <- parent$tkvar
	if(is.null(parent.frm))
		this$tkvar <- tktoplevel()
	else {
		args <- c(list("labelframe", parent.frm), this$conf)
		this$tkvar <- do.call(tkwidget, args)
	}
  nds <- this$childs
  if(!is.null(nds))
  	for(nd in nds)
    	build(nd)
	if(!is.null(parent.frm)) {
  	cmd <- list(parent$cmd, this$tkvar)
  	opts <- parent$opts
  	do.call(tcl, c(cmd, opts))
	}
}


#Creates all necessary (and properly linked) tcltk widgets, and shows them
#'Widget' method
build.Widget <- function(this, ...) {
	parent.frm <- this$parent$tkvar
	args <- c(list(parent.frm), this$opts)
  tmp <- do.call(tkwidget, args)
  this$tkvar <- tmp
  cmd <- list(this$parent$cmd, this$tkvar)
  opts <- this$parent$opts
  do.call(tcl, c(cmd, opts))
}

#From a named list of actions, instantiates a Frame with a button for each listed action
makeButtonsFrame <- function(actions) {
	names <- names(actions)
	fr <- Frame(opts=list(side="left"))
	for(i in 1:length(actions))
		add(fr,Widget(opts=list(type="button", text=names[i], command=actions[[i]])))
	return(fr)
}

#From a root frame, *builds* a complete dialog and shows it
buildDialog <- function(title, rootNode) {
	build(rootNode)
	tkwm.title(rootNode$tkvar, title)
}