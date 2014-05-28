/**
\mainpage libTarget documentation
\section Introduction

The goal of this library is to find targets in images.

It can deal with three kinds of targets:

<table>
<tr>
  <th>OneBlob</th>
  <th>ThreeBlobs</th>
  <th>TwoRings</th>
</tr>
<tr>
  <td>\image html OneBlob.png</td>
  <td>\image html ThreeBlobs.png</td>
  <td>\image html TwoRings.png</td>
</tr>
</table>

Each target include a header, a message and an optional parity bit as follow:

<table>
<tr>
<td>\image html bits.png</td>
<td>
Legend:
<ul> 
<li>Green: header </li> 
<li>Red: message </li> 
<li>Blue: parity (optinal) </li> 
</ul> 
</td>
</tr>
</table>

The target detector is looking for a particular header value given by the user and returns all targets found and their message.

*/
