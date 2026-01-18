
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p> This is the homepage of the <b>Structure Optimized Proximity Scaling (STOPS)</b> project that feature various extended and flexible Multidimensional Scaling approaches. On this page you can find links to related papers, talks, data and software.</p>


<p><a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />Except noted otherwise, content on this homepage and this work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.</p>


<h3>Papers:</h3>
<h4>Peer-reviewed Articles:</h4>
(all open access)</br>
<p>
<dl>
<li>Rusch, T. (2025). <a href="https://onlinelibrary.wiley.com/doi/10.1002/sta4.70086">Multidimensional Scaling With Heaviside Weighting: Extensions to Curvilinear Component Analysis and Curvilinear Distance Analysis</a>Stat, 14:3,e70086.</li>
<li>Rusch, T., Mair, P., & Hornik, K. (2023). <a href="https://link.springer.com/article/10.1007/s11222-022-10197-w">Structure-based Hyperparameter Selection with Bayesian optimization in Multidimensional Scaling</a>, Statistics and Computing, 33:28.</li>
<li>Rusch, T., Mair, P., & Hornik, K. (2021) <a href="https://www.tandfonline.com/doi/full/10.1080/10618600.2020.1869027">Cluster Optimized Proximity Scaling</a>, Journal of Computational and Graphical Statistics, 30:4, 1156-1167.</li>
<li>Rusch, T., Hornik, K. & Mair, P. (2018) <a href="https://www.tandfonline.com/doi/full/10.1080/10618600.2017.1349664">Assessing and Quantifying Clusteredness: The OPTICS Cordillera</a>, Journal of Computational and Graphical Statistics, 27:1, 220-233.</li>
</dl>
</p>

<h4>Technical Reports:</h4>
<p>
<a href="https://research.wu.ac.at/en/publications/sparsified-multidimensional-scaling-and-sparsified-multidimension/"">Sparsified Multidimensional Scaling and Sparsified Multidimensional Distance Analysis</a></br> 
<a href="https://research.wu.ac.at/en/publications/stops-structure-optimized-proximity-scaling-4">STOPS: Structure Optimized Proximity Scaling</a></br> 
<a href="http://epub.wu.ac.at/4888/">COPS: Cluster Optimized Proximity Scaling</a></br> 
<a href="http://epub.wu.ac.at/4789/">Assessing and Quantifying Clusteredness: The OPTICS Cordillera</a> 
</p>

<h3>Talks:</h3>
(title links lead to slides)</br>
<table>
<tr>
  <th>Title</th>
  <th>Event</th>
  <th>Date</th>
  <th>Place</th>
</tr> 
<tr>
  <td><a href="Rusch-ViennaVis25.pdf">Manifold Learning with Extensions to CLCA and CLDA</a></td>
  <td>Vienna Workshop of Model and Data Visualisation</td>
  <td>01.07.2025</td>
  <td>TU Vienna, Austria</td>
</tr> 
<tr>
  <td><a href="stops-psychoco17.pdf">Structure Optimized Proximity Scaling (STOPS): A Framework for Hyperparameter Selection in Multidimensional Scaling</a></td>
  <td>Psychoco 2017</td>
  <td>09.02.2017-10.02.2017</td>
  <td>WU Vienna, Austria</td>
</tr> 
<tr>
  <td><a href="stopsBBS.pdf">COPS and STOPS: Cluster and/or Structure Optimized Proximity Scaling</a></td> 
  <td>Brown Bag Seminar, Institute for Statistics and Mathematics </td>
  <td>07.12.2016</td>
  <td>WU Vienna, Austria</td>
</tr> 
<tr>
  <td><a href="cordilleraBBS.pdf">The OPTICS Cordillera: Nonparametric Assessment of Clusteredness</a></td>
  <td>Brown Bag Seminar, Institute for Statistics and Mathematics </td>
  <td>23.10.2016</td>
  <td>WU Vienna, Austria</td>
</tr> 
<tr>
  <td><a href="http://epub.wu.ac.at/4478/">COPS: Cluster Optimized Proximity Scaling</a></td>
  <td>Psychoco 2015</td>
  <td>12.02.2015-13.02.2015</td>
  <td>Amsterdam, The Netherlands</td>
</tr> 
<tr>
  <td><a href="http://epub.wu.ac.at/4477/">Scaling for Clusters with COPS: Cluster Optimized Proximity Scaling</a>
  <td>CFE-ERCIM 2014</td>
  <td>06.12.2014-08.12.2014</td>
  <td>Pisa, Italy</td>
</tr> 
</table>

<h3>Software:</h3>
<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

The most recent builds of the software packages is available for Windows and Linux here: <a href="https://r-forge.r-project.org/R/?group_id=2037">STOPS Package</a>


<h3>People:</h3>
<dl>
<li><a href="http://www.wu.ac.at/methods/team/dr-thomas-rusch/en/">Thomas Rusch</a></li> 
<li><a href="http://gifi.stat.ucla.edu/">Jan de Leeuw</a></li>
<li><a href="http://scholar.harvard.edu/mair/home">Patrick Mair</a></li>
<li><a href="http://www.wu.ac.at/statmath/en/faculty_staff/faculty/khornik">Kurt Hornik</a></li> 
</dl>
</body>
</html>


