	ۈ'�Y
@ۈ'�Y
@!ۈ'�Y
@	p�Ն @p�Ն @!p�Ն @"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$ۈ'�Y
@����>�?A�#��H	@Y�'�.��?*	P��n]@2j
3Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat�z�΅�?!�4�;�}=@)��~m���?1��io�:@:Preprocessing2t
=Iterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate�;�(A�?!�z�r=@)�E�J�?1�2�wƥ8@:Preprocessing2F
Iterator::Model�7�{�5�?!K�ȷ�>@)Q�O�Iҕ?13e��\2@:Preprocessing2S
Iterator::Model::ParallelMap���%�2�?!0��Z$�(@)���%�2�?10��Z$�(@:Preprocessing2X
!Iterator::Model::ParallelMap::ZipN
�g��?!���~VQ@)hx�輸?1Ƒf�l@:Preprocessing2�
MIterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice	q���v?!� �k�3@)	q���v?1� �k�3@:Preprocessing2v
?Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat::FromTensor�ꫫ�h?!������@)�ꫫ�h?1������@:Preprocessing2d
-Iterator::Model::ParallelMap::Zip[0]::FlatMap�3���?!D��� @@)�y�'Lh?10�Ų[r@:Preprocessing:�
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
�Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
�Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
�Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
�Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)�
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis�
device�Your program is NOT input-bound because only 2.0% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no#You may skip the rest of this page.B�
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown�
	����>�?����>�?!����>�?      ��!       "      ��!       *      ��!       2	�#��H	@�#��H	@!�#��H	@:      ��!       B      ��!       J	�'�.��?�'�.��?!�'�.��?R      ��!       Z	�'�.��?�'�.��?!�'�.��?JCPU_ONLY