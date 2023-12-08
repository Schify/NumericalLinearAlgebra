% Copyright (c) 2011, Dirk-Jan Kroon
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function N=LineNormals2D(Vertices,Lines)
% This function calculates the normals, of the line points
% using the neighbouring points of each contour point, and 
% forward an backward differences on the end points
%
% N=LineNormals2D(V,L)
%
% inputs,
%   V : List of points/vertices 2 x M
% (optional)
%   Lines : A N x 2 list of line pieces, by indices of the vertices
%         (if not set assume Lines=[1 2; 2 3 ; ... ; M-1 M])
%
% outputs,
%   N : The normals of the Vertices 2 x M
%
% Example, Hand
%  load('testdata');
%  N=LineNormals2D(Vertices,Lines);
%  figure,
%  plot([Vertices(:,1) Vertices(:,1)+10*N(:,1)]',[Vertices(:,2) Vertices(:,2)+10*N(:,2)]');
%
% Function is written by D.Kroon University of Twente (August 2011)

% If no line-indices, assume a x(1) connected with x(2), x(3) with x(4) ...
if(nargin<2)
    Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
end

% Calculate tangent vectors
DT=Vertices(Lines(:,1),:)-Vertices(Lines(:,2),:);

% Make influence of tangent vector 1/Distance
% (Weighted Central Differences. Points which are closer give a 
% more accurate estimate of the normal)
LL=sqrt(DT(:,1).^2+DT(:,2).^2);
DT(:,1)=DT(:,1)./max(LL.^2,eps);
DT(:,2)=DT(:,2)./max(LL.^2,eps);

D1=zeros(size(Vertices)); D1(Lines(:,1),:)=DT;
D2=zeros(size(Vertices)); D2(Lines(:,2),:)=DT;
D=D1+D2;

% Normalize the normal
LL=sqrt(D(:,1).^2+D(:,2).^2);
N(:,1)=-D(:,2)./LL;
N(:,2)= D(:,1)./LL;
