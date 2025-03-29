using UnityEngine;

[RequireComponent(typeof(Camera))]
public class FluidRenderEffect : MonoBehaviour
{
    [SerializeField] private Material _compositeMaterial;
    
    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        if (_compositeMaterial != null)
        {
            Graphics.Blit(source, destination, _compositeMaterial);
        }
        else
        {
            Graphics.Blit(source, destination);
        }
    }
}